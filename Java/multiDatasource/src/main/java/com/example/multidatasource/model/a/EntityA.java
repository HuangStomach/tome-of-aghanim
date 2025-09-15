
package com.example.multidatasource.model.a;

import lombok.Data;
import javax.persistence.*;

@Data
@Entity
@Table(name = "entity_a")
public class EntityA {
    @Id
    @GeneratedValue(strategy = GenerationType.IDENTITY)
    private Long id;
    
    private String name;
    private String description;
}
