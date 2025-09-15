
package com.example.multidatasource.repository.b;

import com.example.multidatasource.model.b.EntityB;
import org.springframework.data.jpa.repository.JpaRepository;
import org.springframework.stereotype.Repository;

@Repository
public interface EntityBRepository extends JpaRepository<EntityB, Long> {
}
